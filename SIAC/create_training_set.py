import scipy.stats as stats
from SIAC.lhd import lhd
def create_training_set ( parameters, minvals, maxvals, n_train=200 ):
    """Creates a traning set for a set of parameters specified by 
    ``parameters`` (not actually used, but useful for debugging
    maybe). Parameters are assumed to be uniformly distributed
    between ``minvals`` and ``maxvals``. ``n_train`` input parameter
    sets will be produced, and returned with the actual distributions
    list. The latter is useful to create validation sets.
    Parameters
    -------------
    parameters: list
        A list of parameter names
    minvals: list
        The minimum value of the parameters. Same order as ``parameters``
    maxvals: list
        The maximum value of the parameters. Same order as ``parameters``
    n_train: int
        How many training points to produce
    Returns
    ---------
    The training set and a distributions object that can be used by
    ``create_validation_set``-- Jose:
    https://github.com/jgomezdans/gp_emulator/blob/master/gp_emulator/emulation_helpers.py
    """

    distributions = []
    for i,p in enumerate(parameters):
        distributions.append ( stats.uniform ( loc=minvals[i], \
                            scale=(maxvals[i]-minvals[i] ) ) )
    samples = lhd ( dist=distributions, size=n_train )
    return samples, distributions
