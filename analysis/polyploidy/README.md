extract DAPI only channel from .tif files

polyploidy_extractDAPIonlystacks.txt

	----Quality Control checks-----

	check probability maps generated in ilastik

	polyploidy_QCprobabilitymapcheck.txt

	check alignment of DAPI .tifs and probability .h5 

	polyploidy_QCcheckDAPIandH5alignment.txt

	check segmentation and adjust threshold 

	polyploidy_QCtestnucleisegmentationthreshold.txt

	check intensity measurement and adjust small object filter

	polyploidy_QCDAPIintensitytestsmallobjectfilteredsinglefile.txt

---Quantify images---

run DAPI measurements on bulk files

polyploidy_bulkmeasureDAPI.txt 

collate per image csvs into summary csv

polyploidy_collateintosummarytable.txt

---Stats---

Run in Rstudio

polyploidy_statisticalanalysisR.txt

