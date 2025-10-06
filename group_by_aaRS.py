import pandas as pd

# Load the existing CSV file into a DataFrame
df = pd.read_csv('/Users/chintansoni/Desktop/NGS/TWIST analysis/Rb2NH2/Results/ncAA out>0/RB2_NH2_avg_en_out>0_PB_-ncAA_at2min2.csv')

#grouped_data = df.groupby('sites').agg({'barcode': list, 'average_noncAA': list}).reset_index()

grouped_data = df.groupby('aaRS').agg({'barcode': list}).reset_index()


# Group the data by 'aaRS' and count occurrences
aaRS_counts = df['aaRS'].value_counts().reset_index()
print (aaRS_counts)

# Rename the columns for clarity
aaRS_counts.columns = ['aaRS', 'count']

# Save the grouped data to a new CSV file
grouped_data.to_csv('/Users/chintansoni/Desktop/NGS/TWIST analysis/Rb2NH2/Results/ncAA out>0/RB2_NH2_avg_en_out>0_PB_-ncAA_at2min2_grouped.csv', index=False)
aaRS_counts.to_csv('/Users/chintansoni/Desktop/NGS/TWIST analysis/Rb2NH2/Results/ncAA out>0/RB2_NH2_avg_en_out>0_PB_-ncAA_at2min2_avgd.csv', index=False)