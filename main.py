import pandas as pd


def load_data_from_csv(csv_file):
    print('\nLoading dataset %s' % csv_file.split('/')[-1])
    data = pd.read_csv(csv_file)

    print('shape:', data.shape)


if __name__ == '__main__':
    DATA_DIR = './data/'

    load_data_from_csv(DATA_DIR + 'E13_het.csv')
    load_data_from_csv(DATA_DIR + 'E13_hom.csv')
    load_data_from_csv(DATA_DIR + 'E14_het.csv')
    load_data_from_csv(DATA_DIR + 'E14_hom.csv')
