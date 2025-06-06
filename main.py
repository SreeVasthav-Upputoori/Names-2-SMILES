import requests
import pandas as pd
import concurrent.futures
import time
import random
import logging
import json
import os
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import argparse
from dataclasses import dataclass
from urllib.parse import quote_plus
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('name2smiles.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class Config:
    """Configuration class for the SMILES converter"""
    max_workers: int = 10
    batch_size: int = 5
    timeout: int = 15
    max_retries: int = 5
    rate_limit_delay: float = 0.5
    exponential_backoff_base: float = 2.0
    cache_file: str = "smiles_cache.json"
    proxies: Dict[str, Optional[str]] = None
    
    def __post_init__(self):
        if self.proxies is None:
            self.proxies = {"http": None, "https": None}

class SMILESConverter:
    """Enhanced SMILES converter with caching, rate limiting, and multiple APIs"""
    
    def __init__(self, config: Config = None):
        self.config = config or Config()
        self.cache = self._load_cache()
        self.stats = {
            'total_processed': 0,
            'cache_hits': 0,
            'pubchem_success': 0,
            'opsin_success': 0,
            'failures': 0
        }
    
    def _load_cache(self) -> Dict[str, str]:
        """Load cached SMILES data"""
        try:
            if os.path.exists(self.config.cache_file):
                with open(self.config.cache_file, 'r') as f:
                    cache = json.load(f)
                logger.info(f"Loaded {len(cache)} entries from cache")
                return cache
        except Exception as e:
            logger.warning(f"Failed to load cache: {e}")
        return {}
    
    def _save_cache(self):
        """Save cache to file"""
        try:
            with open(self.config.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
            logger.info(f"Saved {len(self.cache)} entries to cache")
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")
    
    def _clean_name(self, name: str) -> str:
        """Clean and normalize chemical names"""
        if not isinstance(name, str):
            return str(name)
        
        # Remove extra whitespace and normalize
        name = re.sub(r'\s+', ' ', name.strip())
        
        # Handle multiple components (take first)
        name = re.split(r'[;,]', name)[0].strip()
        
        # Remove common prefixes/suffixes that might interfere
        prefixes_to_remove = ['drug:', 'compound:', 'chemical:']
        for prefix in prefixes_to_remove:
            if name.lower().startswith(prefix):
                name = name[len(prefix):].strip()
        
        return name
    
    def _rate_limit(self):
        """Apply rate limiting"""
        time.sleep(self.config.rate_limit_delay + random.uniform(0, 0.3))
    
    def _make_request(self, url: str, retries: int = None) -> Optional[requests.Response]:
        """Make HTTP request with retry logic"""
        retries = retries or self.config.max_retries
        
        for attempt in range(retries):
            try:
                response = requests.get(
                    url, 
                    proxies=self.config.proxies, 
                    timeout=self.config.timeout,
                    headers={'User-Agent': 'SMILES-Converter/1.0'}
                )
                if response.status_code == 200:
                    return response
                elif response.status_code == 404:
                    logger.debug(f"Resource not found: {url}")
                    return None
                elif response.status_code == 429:
                    wait_time = 10 * (self.config.exponential_backoff_base ** attempt)
                    logger.warning(f"Rate limited, waiting {wait_time}s")
                    time.sleep(wait_time)
                else:
                    logger.warning(f"HTTP {response.status_code} for {url}")
                    
            except requests.exceptions.Timeout:
                logger.warning(f"Timeout for {url} (attempt {attempt + 1})")
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request error for {url}: {e} (attempt {attempt + 1})")
            
            if attempt < retries - 1:
                wait_time = self.config.rate_limit_delay * (self.config.exponential_backoff_base ** attempt)
                time.sleep(wait_time)
        
        return None
    
    def _get_smiles_pubchem_batch(self, names: List[str]) -> Dict[str, str]:
        """Get SMILES from PubChem using batch requests"""
        results = {}
        
        try:
            # URL encode names properly
            names_query = ",".join([quote_plus(name) for name in names])
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{names_query}/property/IsomericSMILES/JSON"
            
            response = self._make_request(url)
            if response:
                data = response.json()
                if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                    # Map results back to original names
                    for i, prop in enumerate(data["PropertyTable"]["Properties"]):
                        if i < len(names) and "IsomericSMILES" in prop:
                            results[names[i]] = prop["IsomericSMILES"]
                            
        except Exception as e:
            logger.error(f"PubChem batch error: {e}")
        
        return results
    
    def _get_smiles_pubchem_single(self, name: str) -> Optional[str]:
        """Get SMILES from PubChem for a single compound"""
        try:
            encoded_name = quote_plus(name)
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/IsomericSMILES/JSON"
            
            response = self._make_request(url)
            if response:
                data = response.json()
                if ("PropertyTable" in data and 
                    "Properties" in data["PropertyTable"] and 
                    len(data["PropertyTable"]["Properties"]) > 0):
                    return data["PropertyTable"]["Properties"][0].get("IsomericSMILES")
                    
        except Exception as e:
            logger.debug(f"PubChem single lookup error for '{name}': {e}")
        
        return None
    
    def _get_smiles_opsin(self, name: str) -> Optional[str]:
        """Get SMILES from OPSIN service"""
        try:
            encoded_name = quote_plus(name)
            url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.json"
            
            response = self._make_request(url)
            if response:
                data = response.json()
                if "smiles" in data and data["smiles"]:
                    return data["smiles"]
                    
        except Exception as e:
            logger.debug(f"OPSIN error for '{name}': {e}")
        
        return None
    
    def _get_smiles_cactus(self, name: str) -> Optional[str]:
        """Get SMILES from NCI CACTUS service"""
        try:
            encoded_name = quote_plus(name)
            url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded_name}/smiles"
            
            response = self._make_request(url)
            if response:
                smiles = response.text.strip()
                if smiles and not smiles.lower().startswith('error'):
                    return smiles
                    
        except Exception as e:
            logger.debug(f"CACTUS error for '{name}': {e}")
        
        return None
    
    def get_smiles(self, name: str) -> Optional[str]:
        """Get SMILES for a single compound using multiple sources"""
        if not name or pd.isna(name):
            return None
        
        cleaned_name = self._clean_name(name)
        if not cleaned_name:
            return None
        
        # Check cache first
        cache_key = cleaned_name.lower()
        if cache_key in self.cache:
            self.stats['cache_hits'] += 1
            return self.cache[cache_key]
        
        self.stats['total_processed'] += 1
        smiles = None
        
        # Try PubChem first
        self._rate_limit()
        smiles = self._get_smiles_pubchem_single(cleaned_name)
        if smiles:
            self.stats['pubchem_success'] += 1
            logger.debug(f"PubChem success for '{cleaned_name}'")
        else:
            # Try OPSIN
            self._rate_limit()
            smiles = self._get_smiles_opsin(cleaned_name)
            if smiles:
                self.stats['opsin_success'] += 1
                logger.debug(f"OPSIN success for '{cleaned_name}'")
            else:
                # Try CACTUS as last resort
                self._rate_limit()
                smiles = self._get_smiles_cactus(cleaned_name)
                if smiles:
                    logger.debug(f"CACTUS success for '{cleaned_name}'")
        
        # Cache result (even if None to avoid repeated failures)
        self.cache[cache_key] = smiles
        
        if not smiles:
            self.stats['failures'] += 1
            logger.debug(f"Failed to find SMILES for '{cleaned_name}'")
        
        return smiles
    
    def process_batch(self, names: List[str]) -> Dict[str, Optional[str]]:
        """Process a batch of chemical names"""
        results = {}
        
        # Filter out cached entries
        uncached_names = []
        for name in names:
            cleaned_name = self._clean_name(name)
            cache_key = cleaned_name.lower()
            if cache_key in self.cache:
                results[name] = self.cache[cache_key]
                self.stats['cache_hits'] += 1
            else:
                uncached_names.append(name)
        
        if not uncached_names:
            return results
        
        # Process uncached names with threading
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            future_to_name = {executor.submit(self.get_smiles, name): name for name in uncached_names}
            
            for future in concurrent.futures.as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    results[name] = future.result()
                except Exception as e:
                    logger.error(f"Error processing '{name}': {e}")
                    results[name] = None
        
        return results
    
    def process_dataframe(self, df: pd.DataFrame, name_column: str = 'Name') -> pd.DataFrame:
        """Process a DataFrame with chemical names"""
        if name_column not in df.columns:
            raise ValueError(f"Column '{name_column}' not found. Available columns: {list(df.columns)}")
        
        # Get unique names to avoid duplicate processing
        unique_names = df[name_column].dropna().astype(str).unique().tolist()
        logger.info(f"Processing {len(unique_names)} unique chemical names...")
        
        # Process in batches
        all_results = {}
        batch_size = 50  # Larger batches for threading
        
        for i in range(0, len(unique_names), batch_size):
            batch = unique_names[i:i + batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(unique_names)-1)//batch_size + 1}")
            
            batch_results = self.process_batch(batch)
            all_results.update(batch_results)
            
            # Save cache periodically
            if i % (batch_size * 5) == 0:
                self._save_cache()
        
        # Map results back to DataFrame
        df = df.copy()
        df['SMILES'] = df[name_column].map(all_results)
        
        # Final cache save
        self._save_cache()
        
        return df
    
    def get_statistics(self) -> Dict:
        """Get processing statistics"""
        total = self.stats['total_processed'] + self.stats['cache_hits']
        success_rate = 0
        if total > 0:
            success_rate = ((self.stats['pubchem_success'] + self.stats['opsin_success'] + self.stats['cache_hits']) / total) * 100
        
        return {
            **self.stats,
            'success_rate': f"{success_rate:.1f}%",
            'cache_size': len(self.cache)
        }

def validate_smiles(smiles: str) -> bool:
    """Basic SMILES validation"""
    if not smiles or not isinstance(smiles, str):
        return False
    
    # Basic checks
    if len(smiles) < 2:
        return False
    
    # Check for balanced brackets/parentheses
    brackets = {'(': ')', '[': ']'}
    stack = []
    for char in smiles:
        if char in brackets:
            stack.append(char)
        elif char in brackets.values():
            if not stack or brackets.get(stack.pop()) != char:
                return False
    
    return len(stack) == 0

def get_user_input():
    """Get input and output file names from user"""
    print("=== Chemical Name to SMILES Converter ===")
    print()
    
    # Get input file
    while True:
        input_file = input("Enter the input file path (CSV or Excel): ").strip()
        if not input_file:
            print("Please enter a valid file path.")
            continue
            
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"File not found: {input_file}")
            continue
            
        if input_path.suffix.lower() not in ['.csv', '.xlsx', '.xls']:
            print("Please provide a CSV or Excel file (.csv, .xlsx, .xls)")
            continue
            
        break
    
    # Get output file name
    while True:
        print(f"\nInput file: {input_file}")
        output_file = input("Enter the output Excel file name (e.g., results.xlsx): ").strip()
        
        if not output_file:
            print("Please enter a valid output file name.")
            continue
            
        # Ensure .xlsx extension
        if not output_file.lower().endswith('.xlsx'):
            output_file += '.xlsx'
            
        # Check if file exists and ask for confirmation
        if Path(output_file).exists():
            overwrite = input(f"File '{output_file}' already exists. Overwrite? (y/n): ").strip().lower()
            if overwrite not in ['y', 'yes']:
                continue
                
        break
    
    # Get column name
    print(f"\nOutput file: {output_file}")
    name_column = input("Enter the column name containing chemical names (default: 'Name'): ").strip()
    if not name_column:
        name_column = 'Name'
    
    return input_file, output_file, name_column

def process_file(input_file: str, output_file: str, name_column: str = 'Name', config: Config = None):
    """Process a file containing chemical names"""
    input_path = Path(input_file)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Load data
    logger.info(f"Loading data from {input_file}")
    try:
        if input_path.suffix.lower() == '.csv':
            df = pd.read_csv(input_file)
        elif input_path.suffix.lower() in ['.xlsx', '.xls']:
            df = pd.read_excel(input_file)
        else:
            raise ValueError("Unsupported file format. Use CSV or Excel files.")
    except Exception as e:
        raise ValueError(f"Error reading file: {e}")
    
    logger.info(f"Loaded {len(df)} rows with columns: {list(df.columns)}")
    
    # Check if the specified column exists
    if name_column not in df.columns:
        available_cols = list(df.columns)
        raise ValueError(f"Column '{name_column}' not found. Available columns: {available_cols}")
    
    # Show preview of data
    print(f"\nPreview of data from column '{name_column}':")
    sample_names = df[name_column].dropna().head(5).tolist()
    for i, name in enumerate(sample_names, 1):
        print(f"  {i}. {name}")
    
    if len(df[name_column].dropna()) > 5:
        print(f"  ... and {len(df[name_column].dropna()) - 5} more")
    
    # Ask for confirmation
    proceed = input(f"\nProceed with processing {len(df[name_column].dropna())} chemical names? (y/n): ").strip().lower()
    if proceed not in ['y', 'yes']:
        print("Processing cancelled.")
        return
    
    # Process with SMILES converter
    print("\nStarting SMILES conversion...")
    converter = SMILESConverter(config)
    result_df = converter.process_dataframe(df, name_column)
    
    # Add validation column
    result_df['SMILES_Valid'] = result_df['SMILES'].apply(validate_smiles)
    
    # Save results (always as Excel)
    result_df.to_excel(output_file, index=False)
    
    # Print statistics
    stats = converter.get_statistics()
    print("\n" + "="*50)
    print("PROCESSING COMPLETED!")
    print("="*50)
    
    found_smiles = result_df['SMILES'].notna().sum()
    valid_smiles = result_df['SMILES_Valid'].sum()
    total_compounds = len(result_df)
    
    print(f"Results saved to: {output_file}")
    print(f"Total compounds processed: {total_compounds}")
    print(f"SMILES found: {found_smiles}/{total_compounds} ({found_smiles/total_compounds*100:.1f}%)")
    if found_smiles > 0:
        print(f"Valid SMILES: {valid_smiles}/{found_smiles} ({valid_smiles/found_smiles*100:.1f}% of found)")
    
    print(f"\nDetailed Statistics:")
    for key, value in stats.items():
        print(f"  {key.replace('_', ' ').title()}: {value}")
    
    # Show some examples of results
    print(f"\nSample Results:")
    sample_results = result_df[[name_column, 'SMILES', 'SMILES_Valid']].head(3)
    for _, row in sample_results.iterrows():
        status = "✓" if row['SMILES_Valid'] else "✗" if pd.notna(row['SMILES']) else "?"
        print(f"  {status} {row[name_column]}: {row['SMILES'] if pd.notna(row['SMILES']) else 'Not found'}")
    
    print(f"\nOpen '{output_file}' to view all results.")

def main():
    """Main function - Interactive mode or command line"""
    parser = argparse.ArgumentParser(description="Convert chemical names to SMILES notation")
    parser.add_argument("input_file", nargs='?', help="Input CSV or Excel file (optional)")
    parser.add_argument("-o", "--output", help="Output Excel file (optional)")
    parser.add_argument("-c", "--column", default="Name", help="Column name containing chemical names")
    parser.add_argument("--max-workers", type=int, default=10, help="Maximum number of worker threads")
    parser.add_argument("--batch-size", type=int, default=5, help="Batch size for API requests")
    parser.add_argument("--timeout", type=int, default=15, help="Request timeout in seconds")
    parser.add_argument("--cache-file", default="smiles_cache.json", help="Cache file location")
    parser.add_argument("--clear-cache", action="store_true", help="Clear existing cache")
    parser.add_argument("--interactive", "-i", action="store_true", help="Force interactive mode")
    
    args = parser.parse_args()
    
    # Clear cache if requested
    if args.clear_cache and os.path.exists(args.cache_file):
        os.remove(args.cache_file)
        logger.info("Cache cleared")
    
    # Create config
    config = Config(
        max_workers=args.max_workers,
        batch_size=args.batch_size,
        timeout=args.timeout,
        cache_file=args.cache_file
    )
    
    try:
        # Interactive mode if no input file provided or explicitly requested
        if args.interactive or not args.input_file:
            input_file, output_file, name_column = get_user_input()
        else:
            input_file = args.input_file
            output_file = args.output or (Path(args.input_file).stem + "_with_smiles.xlsx")
            name_column = args.column
            
            # Ensure output is Excel format
            if not output_file.lower().endswith('.xlsx'):
                output_file = Path(output_file).stem + '.xlsx'
        
        process_file(input_file, output_file, name_column, config)
        
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.")
    except Exception as e:
        logger.error(f"Error: {e}")
        print(f"\nError: {e}")
        print("Please check your input file and try again.")

if __name__ == "__main__":
    # Interactive mode when run directly without arguments
    if len(os.sys.argv) == 1:
        try:
            input_file, output_file, name_column = get_user_input()
            
            # Create default config for interactive mode
            config = Config(max_workers=10, timeout=15)
            
            process_file(input_file, output_file, name_column, config)
            
        except KeyboardInterrupt:
            print("\nProcess interrupted by user.")
        except Exception as e:
            print(f"\nError: {e}")
            print("Please check your input file and try again.")
    else:
        main()