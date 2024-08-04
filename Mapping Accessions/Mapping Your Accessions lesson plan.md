# In order to explore your data the first step will be to plot the accessions
## To do this there are a few things that are needed
1. Indentify database from where you will get geolocation data
    - Custom Collection
    - Germplasm Collection (e.g. GRIN-GLOBAL, Genesys)
    - Biogeographic Collection (e.g. GBIF)
    - Herbarium Collection (e.g. Bishop Museum, NY Botanical Garden)
2. Idenitfy your species of interest
3. Download all of the data associated with the species 
5. Identify wild accessions
6. Identify landrace accessions
      - Multiple definitions
          - Arbitraty point in time (e.g. 1970, 1946, 1900)
          - Best definition from Khoury et al., 2022
          -  *"A crop variety or population managed by farmers through cultivation, selection and diffusion, which is typically adapted to a local area and to traditional farming systems, has a recognizable identity and geographic origin, and is often genetically heterogeneous"*
          -  Khoury, C. K., Brush, S., Costich, D. E., Curry, H. A., De Haan, S., Engels, J. M., et al. (2022). Crop genetic erosion: understanding and responding to loss of crop diversity. New Phytologist, 233(1), 84-118.
7. Clean data
      - remove all duplicates (e.g. same accession identifier)
      - remove data that do not have a latitute and longitude
      - count the number of accessions that were collected for the from the same location
8. Save final file as a .CSV with a minimum of three columns (accession name, latitude and longitude)

