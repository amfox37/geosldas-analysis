#!/usr/bin/env python
# coding: utf-8

# Import the library EUMDAC
import eumdac
import datetime
import shutil


# Insert your personal key and secret into the single quotes

consumer_key = '**********************'
consumer_secret = '**********************'

credentials = (consumer_key, consumer_secret)

token = eumdac.AccessToken(credentials)

print(f"This token '{token}' expires {token.expiration}")


# Retrieve all collection objects from DataStore
datastore = eumdac.DataStore(token)

 # Select our collection
selected_collection = datastore.get_collection('EO:EUM:DAT:METOP:SOMO25')


# Set sensing start and end time
start = datetime.datetime(2021, 11, 10, 8, 0)
end = datetime.datetime(2021, 11, 10, 8, 10)

# Retrieve datasets that match our filter
products = selected_collection.search(
    dtstart=start, 
    dtend=end)

print(f'Found Datasets: {len(products)} datasets for the given time range')

for product in products:
    print(str(product))    


for product in products:
    with product.open() as fsrc, \
            open(fsrc.name, mode='wb') as fdst:
        shutil.copyfileobj(fsrc, fdst)
        print(f'Download of product {product} finished.')
print('All downloads are finished.')




