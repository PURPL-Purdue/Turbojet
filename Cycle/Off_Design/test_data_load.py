import pandas as pd

# Load the data and inspect it - FIXED VERSION
map_data = pd.read_csv(
    "reference/etas_import_scaled_full.txt",
    sep=r'\s+',  # Split on any whitespace
    header=None,
    names=['rpm_thousands', 'pressure_ratio', 'efficiency'],
    engine='python'
)

print("Data shape:", map_data.shape)
print("\nFirst 10 rows:")
print(map_data.head(10))
print("\nData types:")
print(map_data.dtypes)
print("\nColumn info:")
print(map_data.info())
print("\nAny NaN values?")
print(map_data.isna().sum())
print("\nPressure ratio stats:")
print(map_data['pressure_ratio'].describe())
print("\nUnique efficiency values:")
print(sorted(map_data['efficiency'].unique()))
