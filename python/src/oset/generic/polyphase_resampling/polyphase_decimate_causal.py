import numpy as np

def polyphase_decimate_causal(in_sig, dec_filter_h, dec_rate, init_cond=None):
    """
    A causal polyphase FIR decimator with continuous block boundaries.

    Parameters
    ----------
    in_sig : array-like
        Input signal (1D list or numpy array).
    dec_filter_h : array-like
        FIR decimation filter's impulse response.
    dec_rate : int
        Decimation factor.
    init_cond : array-like, optional
        Previous input block's trailing samples equal to the filter's memory (same length as dec_filter_h - 1).

    Returns
    -------
    in_sig_decimated : np.ndarray
        Decimated output signal.
    final_cond : np.ndarray
        Last dec_filter_h-1 samples to feed as init_cond for the next block (for continuous processing avoiding edge-effects).
        
    Refs:
        - Vaidyanathan, P. P. (1993). Multirate Systems and Filter Banks. Englewood Cliffs, NJ: Prentice Hall.
        - Crochiere, R. E., & Rabiner, L. R. (1983). Multirate digital signal processing. Prentice Hall.
        
    Reza Sameni, 2025
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET    
    
    """

    x_in = np.asarray(in_sig, dtype=float)
    h = np.asarray(dec_filter_h, dtype=float)
    M = int(dec_rate)
    if M <= 0:
        raise ValueError("dec_rate must be a positive integer")
    L = h.size
    if L == 0:
        return np.array([]), np.array([])

    # State handling
    if init_cond is None:
        z = np.zeros(max(L - 1, 0), dtype=float)
    else:
        z = np.asarray(init_cond, dtype=float)
        if z.size != max(L - 1, 0):
            raise ValueError("init_cond must have length len(dec_filter_h)-1")

    # Concatenate state and block
    x = np.concatenate([z, x_in])
    N = x.size

    # How many valid decimated outputs?
    num_outputs = (N - L) // M + 1
    if num_outputs <= 0:
        return np.array([]), x[-(L - 1):].copy() if L > 1 else np.array([])

    # Polyphase branches: e_m[k] = h[kM + m]
    phases = [h[m::M] for m in range(M)]

    y = np.zeros(num_outputs, dtype=float)
    offset = L - 1  # first output time index in concatenated buffer

    for n in range(num_outputs):
        idx0 = offset + n * M  # "time n x M" in x
        acc = 0.0
        for m in range(M):
            taps = phases[m]
            K = taps.size
            if K == 0:
                continue
            # Indices: idx0 - m - M*np.arange(K)
            idxs = idx0 - m - M * np.arange(K)
            # All idxs are guaranteed in-bounds because of the prepended state
            seg = x[idxs]
            acc += np.dot(taps, seg)
        y[n] = acc

    final_cond = x[-(L - 1):].copy() if L > 1 else np.array([])
    return y, final_cond

def polyphase_decimate_causal_optimized_polyphase(in_sig, dec_filter_h, dec_rate, init_cond=None):
    """
    A causal polyphase FIR decimator with continuous block boundaries, 
    optimized for speed by fully vectorizing the polyphase structure.

    Parameters
    ----------
    in_sig : array-like
        Input signal (1D list or numpy array).
    dec_filter_h : array-like
        FIR decimation filter's impulse response.
    dec_rate : int
        Decimation factor (M).
    init_cond : array-like, optional
        Previous input block's trailing samples (filter memory).

    Returns
    -------
    in_sig_decimated : np.ndarray
        Decimated output signal.
    final_cond : np.ndarray
        Last dec_filter_h-1 samples for the next block.
    """

    x_in = np.asarray(in_sig, dtype=float)
    h = np.asarray(dec_filter_h, dtype=float)
    M = int(dec_rate)
    
    if M <= 0:
        raise ValueError("dec_rate must be a positive integer")
        
    L = h.size # Filter length
    mem_len = max(L - 1, 0) # Filter memory length
    
    if L == 0 or x_in.size == 0:
        return np.array([]), np.array([]) if mem_len > 0 else np.array([])

    # State handling: Prepend memory (z)
    if mem_len == 0:
        z = np.array([])
    elif init_cond is None:
        z = np.zeros(mem_len, dtype=float)
    else:
        z = np.asarray(init_cond, dtype=float)
        if z.size != mem_len:
            raise ValueError(f"init_cond must have length {mem_len} (len(dec_filter_h)-1)")

    # Concatenate state and current block
    x = np.concatenate([z, x_in])
    N = x.size

    # The first index in x corresponding to the first fully-filtered, causal output.
    start_idx = L - 1
    
    # Calculate the total number of decimated outputs that can be produced.
    # The last output index will be at N - 1. We must be safe.
    num_total_outputs = (N - L) // M + 1 
    
    if num_total_outputs <= 0:
        final_cond = x[-mem_len:].copy() if mem_len > 0 else np.array([])
        return np.array([]), final_cond

    # --- Polyphase Vectorization Core ---

    # 1. Prepare Polyphase Filters
    phases = [h[m::M] for m in range(M)]

    # Allocate the final decimated output vector
    y = np.zeros(num_total_outputs, dtype=float)
    
    # Iterate through each of the M polyphase branches
    for m in range(M):
        taps = phases[m]
        K = taps.size # Length of the phase filter
        if K == 0:
            continue
            
        # 2. Determine Required Input Length and Indices
        
        # Calculate the number of input samples needed for this phase to produce 
        # exactly 'num_total_outputs' with a 'valid' convolution of length K.
        len_required = num_total_outputs + K - 1
        
        # Determine the starting index for this phase's input slice: x[L-1 - m].
        start_slice_idx = start_idx - m

        # Construct the full explicit indices for the input segment.
        indices = start_slice_idx + M * np.arange(len_required)
        
        # --- FIX: Robust Index Check ---
        # If the final index exceeds the buffer size (N-1), we must truncate the required length
        # and re-generate the index array to stay within bounds.
        max_valid_idx = N - 1
        
        if indices[-1] > max_valid_idx:
            # Recalculate the true maximum number of decimated outputs possible for this phase.
            # (N - 1 - start_slice_idx) is the number of steps available.
            len_safe_indices = (max_valid_idx - start_slice_idx) // M + 1
            
            # The corresponding length of the convolution input for *this* phase must be:
            len_safe_required = len_safe_indices
            
            # Since the output size 'y' is fixed at num_total_outputs, we must ensure 
            # phase_in is correct for the smallest K. 
            
            # Recalculate len_required based on the safe number of decimated outputs
            num_safe_outputs = len_safe_indices - K + 1
            
            # For the current phase, we use the safe index length *if* it is less than required.
            if len_safe_indices < len_required:
                len_required = len_safe_indices
                indices = start_slice_idx + M * np.arange(len_required)

        # 3. Extract Input Segment
        phase_in = x[indices]
        
        # 4. Perform Vectorized Convolution (Filtering)
        y_m = np.convolve(phase_in, taps, mode='valid')
        
        # 5. Accumulate Results
        # Ensure the output length matches the global num_total_outputs
        if y_m.size > num_total_outputs:
             y_m = y_m[:num_total_outputs]
        elif y_m.size < num_total_outputs:
             # This means the assumption for num_total_outputs was too optimistic 
             # because one phase filter was too long. Truncate the final output y.
             # This is a safe way to handle the last few samples when N/L/M are complex.
             y = y[:y_m.size]
             num_total_outputs = y_m.size
        
        y += y_m

    # The decimated output signal 
    in_sig_decimated = y

    # Final condition: memory for the next block.
    final_cond = x[-mem_len:].copy() if mem_len > 0 else np.array([])
    
    return in_sig_decimated, final_cond