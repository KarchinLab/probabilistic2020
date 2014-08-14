import numpy as np
import logging

logger = logging.getLogger(__name__)

class Bootstrap(object):

    def __init__(self,
                 df,
                 sub_sample,
                 num_iter,
                 col_name='Tumor_Sample'):

        self.df = df
        self.set_sub_sample(sub_sample)
        self.set_num_iter(num_iter)
        self.COLUMN_NAME = col_name
        self.sample_names = self.df[self.COLUMN_NAME].unique()
        self.num_sample_names = len(self.sample_names)
        self.total_count = len(self.df)

    def dataframe_generator(self):
        """Generate subsampled data frames according to the subsample
        and num_samples arguments.
        """

        for i in range(self.num_iter):
            logger.info('Feature generation: Sub-sample rate={0}, Iteration={1} . . .'.format(self.sub_sample, i))
            # get sample names to be used
            prng = np.random.RandomState()
            samp_size = int(self.num_sample_names*self.sub_sample)
            samps_of_interest = prng.choice(self.sample_names, samp_size, replace=True)

            # get data from those sample names
            samp_flag = self.df[self.COLUMN_NAME].apply(lambda x: x in samps_of_interest)
            ixs = samp_flag[samp_flag==True].index
            tmp_df = self.df.ix[ixs].copy()

            logger.info('Finished feature generation: Sub-sample rate={0}, Iteration={1}'.format(self.sub_sample, i))
            yield tmp_df

    def set_subsample(self, subsample):
        """Set the fraction of the original total count to actually sample.

        Sampling is done with replacement. If subsample == 1.0 then it is a
        traditional bootstrap procedure. However, subsamples < 1.0 will
        identify dependencies on database size.

        Parameters
        ----------
        subsample : float
            0 < subsample <= 1.0
        """
        if subsample > 0:
            self.sub_sample = subsample
        else:
            raise ValueError('Subsample should be positive.')

    def set_num_iter(self, num_iter):
        """Set number of iterations for sampling

        Parameters
        ----------
        num_iter : int
            number of iterations for sampling
        """
        if num_iter > 0:
            self.num_iter = num_iter
        else:
            raise ValueError('Number of iterations should be positive.')
