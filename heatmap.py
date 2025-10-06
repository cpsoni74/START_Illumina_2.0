#Wizard: Chintan Soni
#Generating heatmap to compare aaRS activity
#Use heatmap.csv to enter values for fraction of acylation at time t, calculated using Image J analysis of tRNA extension gels.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your CSV data
# The CSV should be in a format where rows and columns represent the heatmap matrix
df = pd.read_csv('/Users/chintansoni/Desktop/START 2.0/Kinetics final/heatmap.csv', index_col=0)  # Assumes first column contains row labels

# Create the heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(df, annot=False, cmap='flare', fmt = '.2f', cbar_kws={'label': '% acylation'})

# Customize labels and title
plt.title('Polyspecificity', fontsize = 16)
plt.xlabel('Monomers', fontsize = 16)
plt.ylabel('aaRS', fontsize = 16)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)

# Display the plot
plt.tight_layout()
plt.show()

#color schemes: flare, rocket, crest