# import modules
from datetime import datetime
from Bio import Entrez
from calendar import monthrange
import pandas as pd
import re
import requests
import json

# Set Entrez email

######################################################################

today = datetime.today()
print(today)