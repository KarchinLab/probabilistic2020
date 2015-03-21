import numpy as np
import logging

logger = logging.getLogger(__name__)

class RandomTumorTypes(object):

    def __init__(self,
                 df,
                 sub_sample,
                 num_iter,
                 table_name='mutation',
                 col_name='Tumor_Type'):
        self.df = df
        self.set_sub_sample(sub_sample)
        self.set_num_iter(num_iter)
        self.TABLE_NAME = table_name
        self.COLUMN_NAME = col_name
        self.tumor_types = self.df[self.COLUMN_NAME].unique()
        self.num_tumor_types = len(self.tumor_types)
        self.total_count = len(self.df)

    def dataframe_generator(self):
        """Generate subsampled data frames according to the sub_sample
        and num_iter arguments.
        """
        #n = int(self.sub_sample * self.total_count)  # number of counts to sample

        for i in range(self.num_iter):
            logger.info('Feature generation: Sub-sample={0}, Iteration={1} . . .'.format(self.sub_sample, i))
            # get sample names to be used
            prng = np.random.RandomState()
            prng.shuffle(self.tumor_types)

            types_of_interest = set(self.tumor_types[:int(self.num_tumor_types*self.sub_sample)])

            # get data from those sample names
            samp_flag = self.df[self.COLUMN_NAME].apply(lambda x: x in types_of_interest)
            ixs = samp_flag[samp_flag==True].index
            tmp_df = self.df.ix[ixs].copy()

            logger.info('Finished feature generation: Sub-sample={0}, Iteration={1}'.format(self.sub_sample, i))
            yield tmp_df

    def set_sub_sample(self, sub_sample):
        """Set the fraction of the original total mutations to actually sample.

        Sampling is done without replacement.

        Parameters
        ----------
        sub_sample : float
            0 < sub_sample <= 1.0
        """
        if 0 <= sub_sample <= 1:
            self.sub_sample = sub_sample
        else:
            raise ValueError('Subsample should be between zero and one.')

    def set_num_iter(self, num_iter):
        """Set number of times to sample w/o replacement.

        Parameters
        ----------
        num_iter : int
            do sample w/o replacement, num_iter number of times
        """
        if iter > 0:
            self.num_iter = num_iter
        else:
            raise ValueError('Number of iterations should be positive.')
