import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import argparse
import os
from statsmodels.stats.multitest import multipletests

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Save plain text from a website.')
parser.add_argument('--url',  help='URL of the website', default = "http://159.149.160.88/pscan/")
parser.add_argument('--file', required=True, help='Path to file containing REFSEQID genes')
parser.add_argument('--headless', help = 'bool if a browser opens', default = 1, type = int)
args = parser.parse_args()

#selenium_service = Service('/usr/local/bin/chromedriver')
#service = Service()

df = pd.read_csv(args.file, sep=',')
df_wo_NA = df.dropna()
file_contents = '\\n'.join(df_wo_NA['refseq_ids'].astype(str))


word_count = len(df_wo_NA)
print("Finished Pscan on {}".format(os.path.basename(args.file)))
print("Number of genes:", word_count)
options = webdriver.ChromeOptions()

if (int(args.headless)):
    options.add_argument('--headless')  # Run Chrome in headless mode

#driver = webdriver.Chrome( options=chrome_options)
# driver = webdriver.Chrome(service=service, options=options)
driver = webdriver.Chrome( options=options)

# Set the URL where the form is located
url = args.url

# Load the page
driver.get(url)

# Find and fill the input fields
id_text_input = driver.find_element(By.NAME, 'id_text')
#id_text_input.send_keys(file_contents)
cmd = "document.getElementById('id_text').value = '{}';".format(file_contents)
driver.execute_script(cmd)

# Submit the form
submit_button = driver.find_element(By.CSS_SELECTOR, 'input[type="submit"]')
submit_button.click()

# Find the download link
download_link = driver.find_element(By.LINK_TEXT, 'View Text Results')
download_link.click()

# Switch to the new window
driver.switch_to.window(driver.window_handles[-1])

# Get the plain text from the new window
plain_text = driver.find_element(By.TAG_NAME, 'pre').text

if not os.path.exists("Pscan_results"):
    os.makedirs("Pscan_results")

outfile = os.path.join("Pscan_results", (os.path.basename(os.path.splitext(args.file)[0]) + "_results.csv"))

# Save the plain text to a file
with open(outfile , 'w', encoding='utf-8') as file:
    file.write(plain_text)

# Close the new window
driver.close()
driver.quit()

results = pd.read_csv(outfile, sep = " ")
results['TF_NAME'] = results['TF_NAME'].str.replace(r'\(var\.[0-9]+\)', '', regex=True)
results['Pscan_TERM'] = results['TF_NAME'].copy()
results.insert(1, 'Pscan_TERM', results.pop('Pscan_TERM'))

results = results.assign(TF_NAME=results.TF_NAME.str.split('::'))
results = results.explode("TF_NAME")

overlap = results['TF_NAME'].str.lower().isin(df['GENE_SYMBOLS'].str.lower())
results['In_original_list'] = overlap
results['FDR'] = multipletests(results['P_VALUE'], method='fdr_bh')[1]

results.to_csv(outfile, index=False)