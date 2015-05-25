% Initialize test dataExtractor

improc2.tests.cleanupForTests;
inMemoryCollection = improc2.tests.data.collectionOfProcessedDAGObjects();

dataExtractor = improc2.launchDataExtractor(inMemoryCollection);

% Run following command, and look at output in Matlab terminal.
% This will tell you all of the methods associated with the extractor.
dataExtractor

% To see what data is available to you (e.g. in cy), run this:
availableData = dataExtractor.getDataFromObjectHandle('cy');

% To take a look at it, just type:
availableData

% Alternatively, type "availableData." and hit the 'Tab' button - a list of
% functions will pop up. Things like "getNumSpots" and "zMerge". Play
% around with these and see what all you have available to you.

% For example, the following command should give you the number of cy spots
% in the first object in your first data file:
availableData.getNumSpots()

% You can also access the x, y, and z positions for each spot using:
[x,y,z] = availableData.getSpotCoordinates()

% Now understand the following command (this shows up on the wiki as well):
% The "@getNumSpots" in this command is actually a function that takes
% "availableData" as an argument. That's how it is able to get you the
% number of spots for each object.
dataExtractor.extractFromProcessorData('cyRnaCount', @getNumSpots, 'cy:Spots');

% Check to see that the above line of code actually gave you the right
% number of spots for objArrayNum = 1 (first data file) and objNum = 1
% (first object number)
dataExtractor.extractAllToCSVFile('test.csv')

% With all this in mind, let's see if we can extract the z difference (max
% z - min z) for all our image objects. Check out the "findZdifference.m"
% function. Understand what it's taking as an argument, and how it's giving
% the correct output. (Is it giving the correct output?)

% Can you write an extract code that will allow you to extract the z
% difference for both cy and alexa?

% Finally, using the "availableData" for your own dataset (not this test
% dataset we've been working with so far), can you figure out how to
% extract the number of nuclei?