'''
Example of WGS84 to EGSA87 conversion and the opposite

Dimitrios Dimopoulos 2021
'''

import math
from egsa87wgs84 import egsa87wgs84
from wgs84egsa87 import wgs84egsa87


option = input("Press 1 for WGS84 to EGSA87 or 2 for EGSA87 to WGS84 :")

if option == "1":
    longitude = input("Please enter longitude in rads:")
    latitude = input("Please enter latitude in rads :")
    print(wgs84egsa87(float(longitude),float(latitude)))
elif option =='2':
    x = input("Please enter x coords in meters:")
    y = input("Please enter y coords in meters:")
    print(egsa87wgs84(float(x),float(y)))
else:
    print("Invalid option")
    
