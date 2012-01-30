def call( ctrl, error_parameters, **keyword_args ):
    """
    arguments should be
      addressing( control , epsilon )
      
    outputs should be
      new_control
    """
    
    # Enforce defaults
    if error_parameters == None:
        error_parameters = default_parameters()
      
    # Perhaps the user was lasy any input the error parameters as a
    # single element rather than a subscriptable list.  Fix.
    try:
        epsilon = error_parameters[0]
    except TypeError:
        epsilon = error_parameters
    
    # Update amplitude of control values.
    ctrl.control = (1.0 + epsilon) * ctrl.control
    
    # Return modified control.
    return ctrl


def default_parameters():
    """
    default parameters
    """
    
    return [0.1]


def repr( error_parameters ):
    """
    Function to display amplitude objects when called on the command line
    """
    
    string = "amplitude error: \n" + \
             "\tepsilon:\t%.2E" %( error_parameters[0] )
    return string
