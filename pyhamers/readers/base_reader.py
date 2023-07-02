import abc

class BaseReader(metaclass=abc.ABCMeta):
    """
    Abstract base class to read data.
    
    Steps in using the data reader to post-process data in many time steps:
    1. Use the constructor BaseReader(data_directory_path) to create the object and initialize it.
    2. Get the full domain size.
    3. Call setSubDomain((lo, hi)) to set the sub-domain to read in.
    4. Call readCoordinates() to get coordinates of the sub-domain.
    For each time step:
        a. Call setStep(step) for each timestep.
        b. Call readData().
        c. Do your post-processing...
    
    To write a concrete class (called MyReaderImplementation, say) that derives from this, implement the following
    abstract methods and in the end of the file add the following code to register the concrete class
    BaseReader.register(MyReaderImplementation).
    """
    
    @abc.abstractmethod
    def dimension(self):
        """
        Return the dimension of the domain.
        """
        return
    
    
    @abc.abstractmethod
    def setStep(self, step):
        """
        Update the metadata from the summary file in the data directory at a new time step.
        """
        return
    
    
    @abc.abstractmethod
    def getStep(self):
        """
        Return the time step that is currently set.
        """
        return
    
    
    step = abc.abstractproperty(getStep, setStep)
    
    
    @abc.abstractproperty
    def domain_size(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """
        return
    
    
    @abc.abstractmethod
    def setSubDomain(self, lo_and_hi):
        """
        Set the sub-domain for reading coordinates and data.
        """
        return
    
    
    @abc.abstractmethod
    def getSubDomain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """
        return
    
    
    sub_domain = abc.abstractproperty(getSubDomain, setSubDomain)
    
    
    @abc.abstractproperty
    def periodic_dimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """
        return
    
    
    @abc.abstractproperty
    def time(self):
        """
        Return the simulation time at current time step.
        """
        return
    
    
    @abc.abstractproperty
    def steps(self):
        """
        Return all of the steps.
        """
        return
    
    
    @abc.abstractproperty
    def data_order(self):
        """
        Return the data order.
        """
        return
    
    
    @abc.abstractmethod
    def readCoordinates(self):
        """
        Get the coordinates of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """
        return
    
    
    @abc.abstractmethod
    def readData(self, var_names, data=None):
        """
        Read the data of several variables in the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """
        return