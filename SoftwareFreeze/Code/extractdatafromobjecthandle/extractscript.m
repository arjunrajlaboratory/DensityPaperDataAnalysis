% extract max z position
dataExtractor = improc2.launchDataExtractor();

dataExtractor.extractFromProcessorData('zDiffAlexa', @findZdifference, 'alexa:Spots');

dataExtractor.extractAllToCSVFile('test.csv')


% extract number of nuclei
extractor.extractFromProcessorData('numberOfNuclei', whatToExtract, 'nuclei')

