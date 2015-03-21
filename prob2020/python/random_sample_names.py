import numpy as np
import logging

logger = logging.getLogger(__name__)

class RandomSampleNames(object):

    def __init__(self,
                 df,
                 sub_sample,
                 num_iter,
                 table_name='mutation',
                 col_name='Tumor_Sample'):
        self.df = df
        self.set_sub_sample(sub_sample)
        self.set_num_iter(num_iter)
        self.COLUMN_NAME = col_name
        self.sample_names = self.df[self.COLUMN_NAME].unique()
        self.num_sample_names = len(self.sample_names)
        self.total_count = len(self.df)

    def dataframe_generator(self):
        """Generate subsampled data frames according to the sub_sample
        and num_iter arguments.
        """

        for i in range(self.num_iter):
            logger.info('Feature generation: Sub-sample rate={0}, Iteration={1} . . .'.format(self.sub_sample, i))
            # get sample names to be used
            prng = np.random.RandomState()
            prng.shuffle(self.sample_names)

            samps_of_interest = set(self.sample_names[:int(self.num_sample_names*self.sub_sample)])

            # get data from those sample names
            samp_flag = self.df[self.COLUMN_NAME].apply(lambda x: x in samps_of_interest)
            ixs = samp_flag[samp_flag==True].index
            tmp_df = self.df.ix[ixs].copy()


            # drop gene name column
            # proc_feat_df = proc_feat_df.drop('gene', axis=1)

            logger.info('Finished feature generation: Sub-sample rate={0}, Iteration={1}'.format(self.sub_sample, i))
            yield tmp_df

    def set_sub_sample(self, sub_sample):
        """Set the fraction of the original total mutations to actually sample.

        Sampling is done without replacement.

        **Parameters**

        sub_sample : float
            0 < sub_sample <= 1.0
        """
        if 0 <= sub_sample <= 1:
            self.sub_sample = sub_sample
        else:
            raise ValueError('Subsample should be between zero and one.')

    def set_num_iter(self, num_iter):
        """Set number of times to sample w/o replacement.

        **Parameters**

        num_iter : int
            do sample w/o replacement, num_iter number of times
        """
        if iter > 0:
            self.num_iter = num_iter
        else:
            raise ValueError('Number of iterations should be positive.')

