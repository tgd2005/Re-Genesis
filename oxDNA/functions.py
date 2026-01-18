import numpy as np

def moving_average_filter(data, window_size=10, epoch_col=None, filter_cols=None, suffix='_filtered'):
    """
    Apply a simple moving average filter to smooth signals in a dataset.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the data to filter
    window_size : int
        Number of samples in the moving average window (default: 10)
    epoch_col : str or None
        Name of the epoch/time column to preserve unfiltered. If None, all columns
        except those in filter_cols will be preserved as-is (default: None)
    filter_cols : str, list of str, or None
        Column name(s) to apply the filter to. If None, applies to all numeric columns
        except epoch_col (default: None)
    suffix : str
        Suffix to add to filtered column names (default: '_filtered')
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with original columns plus filtered columns with suffix added
        
    Examples:
    ---------
    # Filter single column with explicit epoch
    df_filtered = moving_average_filter(df, window_size=50, epoch_col='step', filter_cols='force_mag')
    
    # Filter multiple columns
    df_filtered = moving_average_filter(df, window_size=20, filter_cols=['force_x', 'force_y', 'force_z'])
    
    # Filter all numeric columns automatically
    df_filtered = moving_average_filter(df, window_size=30, epoch_col='time')
    """
    data_copy = data.copy()
    
    # Determine which columns to filter
    if filter_cols is None:
        # Auto-detect numeric columns
        numeric_cols = data_copy.select_dtypes(include=[np.number]).columns.tolist()
        # Exclude epoch column if specified
        if epoch_col is not None and epoch_col in numeric_cols:
            numeric_cols.remove(epoch_col)
        cols_to_filter = numeric_cols
    elif isinstance(filter_cols, str):
        # Single column specified as string
        cols_to_filter = [filter_cols]
    else:
        # List of columns
        cols_to_filter = list(filter_cols)
    
    # Apply rolling mean to each specified column
    for col in cols_to_filter:
        if col in data_copy.columns:
            filtered_col_name = f"{col}{suffix}"
            data_copy[filtered_col_name] = data_copy[col].rolling(
                window=window_size,
                center=True,
                min_periods=1
            ).mean()
        else:
            print(f"Warning: Column '{col}' not found in dataset. Skipping.")
    
    return data_copy