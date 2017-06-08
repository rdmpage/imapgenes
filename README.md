# imapgenes
Data from "An Anthropocene map of genetic diversity"

Miraldo, A., Li, S., Borregaard, M. K., Florez-Rodriguez, A., Gopalakrishnan, S., Rizvanovic, M., … Nogues-Bravo, D. (2016, September 29). An Anthropocene map of genetic diversity. Science. American Association for the Advancement of Science (AAAS). https://doi.org/10.1126/science.aaf4381

Article and supplementary data (Miraldo.SM.pdf Data S2, and Curated_codes Data S1), data seqs_info_tsv from http://macroecology.ku.dk/resources/imapgenes/seqs_info_tsv.zip (see the iMapsGene web site http://macroecology.ku.dk/resources/imapgenes/imapgenes_app_beta/ ).

## Contacts
David Nogués-Bravo dnogues@snm.ku.dk
Carsten Rahbek crahbek@snm.ku.dk
Andreia Miraldo andreia.miraldo@snm.ku.dk

## Examples

The first sequence in the amphibian data file is [DQ523091](https://www.ncbi.nlm.nih.gov/nuccore/DQ523091):
```
/isolate=“EtrivTSValviii9”
/country=“Peru: San Martin, Tarapoto”
```
If we use the geonames API  http://api.geonames.org/searchJSON/searchJSON?username=demo&q=Peru%3A%20San%20Martin%2C%20Tarapoto&maxRows=10 (you will need your own username) we get:

``` 
{
    “totalResultsCount”: 12,
    “geonames”: [
        {
            “adminCode1”: “22”,
            “lng”: “-76.37325”,
            “geonameId”: 6300820,
            “toponymName”: “Tarapoto”,
            “countryId”: “3932488”,
            “fcl”: “S”,
            “population”: 0,
            “countryCode”: “PE”,
            “name”: “Tarapoto”,
            “fclName”: “spot, building, farm”,
            “countryName”: “Peru”,
            “fcodeName”: “airport”,
            “adminName1”: “San Martín”,
            “lat”: “-6.50874”,
            “fcode”: “AIRP”
        },
        .
        .
        .
}
```
This gives the latitude and longitude values (-6.50874, -76.37325) reported in the data dump. Note that it is the airport.

In the paper [Genetic divergence and speciation in lowland and montane peruvian poison frogs](http://dx.doi.org/10.1016/j.ympev.2006.05.005) that published the sequence there is a specimen 
```
E. trivittatusP14	Tarapoto, San Martin, Peru	540 m	S 6.43066′ W 76.29034′
```
There doesn’t seem to be an obvious match between all the names in the table and the sequence in Genbank, but we should be able to extract some records. 




