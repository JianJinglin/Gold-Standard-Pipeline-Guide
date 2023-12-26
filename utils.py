import os
import io
import pandas as pd
import gzip
import mygene

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sparse_lmm import VariableSelection
from statsmodels.stats.multitest import multipletests


def get_relevant_filepaths(cohort_dir):
    files = os.listdir(cohort_dir)
    soft_files = [f for f in files if 'soft' in f.lower()]
    matrix_files = [f for f in files if 'matrix' in f.lower()]
    assert len(soft_files) > 0 and len(matrix_files) > 0
    # If there are multiple SOFT files or matrix files, simply choose the first one. May be replaced by better strategies later.
    soft_file_path = os.path.join(cohort_dir, soft_files[0])
    matrix_file_path = os.path.join(cohort_dir, matrix_files[0])

    return soft_file_path, matrix_file_path

def line_generator(source, source_type):
    """
    Generator that yields lines from a file or a string.

    Parameters:
    - source: File path or string content.
    - source_type: 'file' or 'string'.
    """
    if source_type == 'file':
        with gzip.open(source, 'rt') as f:
            for line in f:
                yield line.strip()
    elif source_type == 'string':
        for line in source.split('\n'):
            yield line.strip()
    else:
        raise ValueError("source_type must be 'file' or 'string'")

def filter_content_by_prefix(source, prefixes, unselect, source_type, return_df=False):
    """
    Filter lines from a file or a string based on specified prefixes.

    Parameters:
    - source: File path or string content to filter.
    - prefixes: List of prefixes to filter by.
    - unselect: Boolean flag to invert selection.
    - return_df: Boolean flag to return a pandas DataFrame.
    - source_type: 'file' or 'string'.
    """
    filtered_lines = []
    prefix_set = set(prefixes)

    # Use generator to get lines
    for line in line_generator(source, source_type):
        matched = any(line.startswith(prefix) for prefix in prefix_set)
        if matched != unselect:
            filtered_lines.append(line)

    filtered_content = '\n'.join(filtered_lines)

    if return_df:
        filtered_content = pd.read_csv(io.StringIO(filtered_content), delimiter='\t', low_memory=False)

    return filtered_content

def get_background_info(file_path, prefixes=['!Sample_geo_accession', '!Series_title', '!Series_summary', '!Series_overall_design', '!Sample_characteristics_ch1']):
    background_info = filter_content_by_prefix(file_path, prefixes, unselect=False, source_type='file', return_df=False)
    return background_info

def get_clinical_data(content, prefixes=['!Sample_geo_accession', '!Sample_characteristics_ch1']):
    clinical_data = filter_content_by_prefix(content, prefixes, unselect=False, source_type='string', return_df=True)
    return clinical_data

def get_gene_annotation(file_path, prefixes=['^', '!', '#']):
    gene_metadata = filter_content_by_prefix(file_path, prefixes, unselect=True, source_type='file', return_df=True)
    return gene_metadata

def get_gene_mapping(annotation, gene_col):
    """ Process the metadata to get the mapping between gene names and gene probes.
    Some of the series only provide accession of gene in GenBank (e.g., GB_ACC), so you need to use GB_ACC to search for the gene name in GeneBank website.
    """
    # mapping_data = filter_lines_by_prefix(file_path, prefixes=['^', '!', '#'], unselect=True, return_df=True)
    # Find the name of the column that stores the gene symbol
    # standardized_cols = mapping_data.columns.str.lower().str.replace(' ', '_')
    #gene_symbol_cols = mapping_data.columns[standardized_cols == 'gene_symbol']
    # assert len(gene_symbol_cols) == 1, f"Expected one 'gene_symbol' column, found {len(gene_symbol_cols)}"
    #gene_symbol_col = gene_symbol_cols[0]
    #mapping_data['ID'] = mapping_data['ID'].astype('int')
    mapping_data = annotation.loc[:, ['ID', gene_col]]
    mapping_data = mapping_data.dropna()
    mapping_data = mapping_data.rename(columns={gene_col: 'Gene'}).astype({'ID': 'str'})

    return mapping_data

def get_genetic_data(file_path):
    genetic_data = pd.read_csv(file_path, compression='gzip', skiprows=52, comment='!', delimiter='\t')
    genetic_data = genetic_data.dropna()
    genetic_data = genetic_data.rename(columns={'ID_REF': 'ID'}).astype({'ID': 'str'})

    return genetic_data

def apply_gene_mapping(expression_df, mapping_df):
    mapping_df['Gene'] = mapping_df['Gene'].str.split(';')
    mapping_df = mapping_df.explode('Gene')
    expression_df = pd.merge(mapping_df, expression_df, on='ID').drop(columns='ID').set_index('Gene')
    expression_df = expression_df.groupby('Gene').mean()

    return expression_df

def normalize_gene_symbols(gene_symbols, batch_size=1000):
    mg = mygene.MyGeneInfo()
    normalized_genes = {}

    # Process in batches
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i + batch_size]
        results = mg.querymany(batch, scopes='symbol', fields='symbol', species='human')

        # Update the normalized_genes dictionary with results from this batch
        for gene in results:
            normalized_genes[gene['query']] = gene.get('symbol', None)

    # Return the normalized symbols in the same order as the input
    return [normalized_genes.get(symbol) for symbol in gene_symbols]

def normalize_gene_symbols_in_index(gene_df):
    normalized_gene_list = normalize_gene_symbols(gene_df.index.tolist())
    assert len(normalized_gene_list) == len(gene_df.index)
    gene_df.index = normalized_gene_list
    gene_df = gene_df[gene_df.index.notnull()]
    return gene_df

def get_feature_data(clinical_df, row_id, feature, decode_fn):
    clinical_df = clinical_df.iloc[row_id:row_id+1, 1:]
    clinical_df.index = [feature]
    clinical_df = clinical_df.applymap(decode_fn)

    return clinical_df

def add_binary_feature(genetic_df, feature_df, feature):
    merged_df = pd.concat([genetic_df, feature_df])
    merged_df = merged_df.loc[:, merged_df.loc[feature].isin([0, 1])]
    return merged_df


def plot_numeric_distribution(df, column):
    plt.figure(figsize=(10, 6))
    sns.histplot(df[column], kde=True, bins=30)
    plt.title(f'Distribution of {column.capitalize()}')
    plt.xlabel('')
    plt.ylabel('Frequency')
    plt.show()


def plot_categorical_distribution(df, column):
    plt.figure(figsize=(10, 6))
    sns.countplot(y=column, data=df, order=df[column].value_counts().index)
    plt.title(f'Distribution of {column.capitalize()}')
    plt.xlabel('Frequency')
    plt.ylabel('')
    plt.show()


def analyze_distributions(df, numerical_columns, categorical_columns):
    for col in numerical_columns:
        plot_numeric_distribution(df, col)

    for col in categorical_columns:
        plot_categorical_distribution(df, col)

def normalize_data(X_train, X_test=None):
    mean = np.mean(X_train, axis=0)
    std = np.std(X_train, axis=0)

    # Handling columns with std = 0
    std_no_zero = np.where(std == 0, 1, std)

    # Normalize X_train
    X_train_normalized = (X_train - mean) / std_no_zero
    # Set normalized values to 0 where std was 0
    X_train_normalized[:, std == 0] = 0

    if X_test is not None:
        X_test_normalized = (X_test - mean) / std_no_zero
        X_test_normalized[:, std == 0] = 0
    else:
        X_test_normalized = None

    return X_train_normalized, X_test_normalized


def cross_validation_with_lasso(X, y, k=5):
    indices = np.arange(X.shape[0])
    np.random.shuffle(indices)

    fold_size = len(X) // k
    accuracies = []

    for i in range(k):
        # Split data into train and test based on the current fold
        test_indices = indices[i * fold_size: (i + 1) * fold_size]
        train_indices = np.setdiff1d(indices, test_indices)

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = y[train_indices], y[test_indices]

        # Preprocess the train and test data
        X_train, X_test = normalize_data(X_train, X_test)

        # Use the LogisticRegression model to train and predict
        model = LogisticRegression(penalty='l1', solver='liblinear', random_state=42)
        model.fit(X_train, y_train)
        predictions = model.predict(X_test)

        # Calculate accuracy for this fold
        accuracy = np.mean(predictions == y_test)
        accuracies.append(accuracy)

    return np.mean(accuracies), np.std(accuracies)

# k-fold cross-validation for the variable selection model
def cross_validation_with_lmm(X, y, k=5):
    indices = np.arange(X.shape[0])
    np.random.shuffle(indices)

    fold_size = len(X) // k
    accuracies = []

    for i in range(k):
        # Split data into train and test based on the current fold
        test_indices = indices[i * fold_size: (i + 1) * fold_size]
        train_indices = np.setdiff1d(indices, test_indices)

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = y[train_indices], y[test_indices]

        # Preprocess the train and test data
        normalized_X_train, normalized_X_test = normalize_data(X_train, X_test)

        # Use the precision_lasso package for regression
        model = VariableSelection()
        model.fit(normalized_X_train, y_train)
        predictions = model.predict(normalized_X_test)

        # Turn the predictions into binary values using a threshold of 0.5
        binary_predictions = (predictions > 0.5).astype(int)

        # Calculate accuracy for this fold
        accuracy = np.mean(binary_predictions == y_test)
        accuracies.append(accuracy)

    return np.mean(accuracies), np.std(accuracies)


def get_feature_related_genes(file_path, feature, normalize):
    related_gene_df = pd.read_csv(file_path)
    related_gene_df = related_gene_df.loc[:, ['Disease Name', 'Corresponding_Gene_Symbol']].set_index('Disease Name')
    feature_related_genes = related_gene_df.loc[feature].tolist()[0].strip().split(',')
    feature_related_genes = [gn.strip() for gn in feature_related_genes]
    if normalize:
        feature_related_genes = normalize_gene_symbols(feature_related_genes)

    return feature_related_genes


def get_gene_regressors(trait_df, condition_df, related_genes):
    gene_regressors = None
    genes_in_trait_data = set(trait_df.index)
    genes_in_condition_data = set(condition_df.index)

    common_genes_across_data = genes_in_trait_data.intersection(genes_in_condition_data)
    if len(common_genes_across_data) == 0:
        print("The trait and condition datasets have no genes in common. Please try other datasets")
    else:
        print(f"The trait and condition datasets have {len(common_genes_across_data)} genes in common, such as {list(common_genes_across_data)[:10]}.")
        common_genes = common_genes_across_data.intersection(set(related_genes))
        if len(common_genes) > 0:
            gene_regressors = list(common_genes)
            print(f"Found {len(common_genes)} candidate genes that can be used in two-step regression analysis, such as {gene_regressors[:10]}.")
        else:
            print(f"The condition and trait datasets have common genes, and we need to select indicator genes for the condition")

    return gene_regressors, common_genes_across_data


def report_result_from_lmm(model, feature_names, trait, condition, threshold=0.05, save_output=True, output_dir='./output'):
    # Retrieve the coefficients
    coefficients = model.getBeta().reshape(-1).tolist()
    nlog_p_values = model.getNegLogP().reshape(-1).tolist()
    p_values = [np.exp(-p) for p in nlog_p_values]

    # Create a DataFrame for the regression results
    regression_df = pd.DataFrame({
        'Variable': feature_names,
        'Coefficient': coefficients,
        'p_value': p_values
    })

    # Extract information about the condition's effect
    condition_effect = regression_df[regression_df['Variable'] == condition].iloc[0]

    # Report the effect of the condition
    print(f"Effect of the condition on the target variable:")
    print(f"Variable: {condition}")
    print(f"Coefficient: {condition_effect['Coefficient']:.4f}")
    print(f"p-value: {condition_effect['p_value']:.4g}\n")

    gene_regression_df = regression_df[regression_df['Variable'] != condition]

    # Apply the Benjamini-Hochberg correction, to get the corrected p-values
    corrected_p_values = multipletests(gene_regression_df['p_value'], alpha=threshold, method='fdr_bh')[1]
    gene_regression_df.loc[:, 'corrected_p_value'] = corrected_p_values

    significant_genes = gene_regression_df.loc[gene_regression_df['corrected_p_value'] < threshold]
    significant_genes_sorted = significant_genes.sort_values('corrected_p_value')

    print(f"Found {len(significant_genes_sorted)} significant genes affecting the trait '{trait}' conditional on the factor '{condition}', with corrected p-value < {threshold}:")

    print(significant_genes_sorted[['Variable', 'Coefficient', 'corrected_p_value']].to_string(index=False))

    # Optionally, save this to a CSV file
    if save_output:
        os.makedirs(output_dir, exist_ok=True)
        significant_genes_sorted.to_csv(os.path.join(output_dir, f'{trait}_{condition}_significant_genes.csv'), index=False)


def judge_binary_variable_biased(dataframe, col_name, min_proportion, min_num):
    label_counter = dataframe[col_name].value_counts()
    total_samples = len(dataframe)
    rare_label_num = label_counter.min()
    rare_label = label_counter.idxmin()
    rare_label_proportion = rare_label_num / total_samples

    print(f"The least common label is '{rare_label}' with {rare_label_num} occurrences. This represents {rare_label_proportion:.2%} of the dataset.")

    biased = (len(label_counter) < 2) or ((rare_label_proportion < min_proportion) and (rare_label_num < min_num))
    return biased


def judge_continuous_variable_biased(dataframe, col_name):
    # Calculate the quartiles of the continuous variable
    quartiles = dataframe[col_name].quantile([0.25, 0.5, 0.75])
    min_value = dataframe[col_name].min()
    max_value = dataframe[col_name].max()

    # Printing quartile information
    print(f"Quartiles for '{col_name}':")
    print(f"  25%: {quartiles[0.25]}")
    print(f"  50% (Median): {quartiles[0.5]}")
    print(f"  75%: {quartiles[0.75]}")
    print(f"Min: {min_value}")
    print(f"Max: {max_value}")

    # Check if the variable is too biased for regression analysis.
    # As a starting point, we consider it biased if all values are the same
    # For the next step, maybe ask GPT to judge based on the quartile statistics combined with its common sense knowledge about this feature.
    biased = min_value == max_value

    return biased