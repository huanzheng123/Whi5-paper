// @File(label = "Input directory", style = "directory") input
a=0
list = getFileList(input);
list = Array.sort(list);


for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], ".TIF")){
    //setTool("line");
    filename = input +"/"+ list[i];
    filenameMod=filename;
    print(filenameMod);

    name=split(filenameMod,".");
    filenameNew= name[0]+"new_disk30.tif";
    open(filenameMod);
    name=split(list[i],".");

	//run("Mean...", "radius=1.0");  //added this denoising
	run("Subtract Background...", "rolling=20 disable stack"); //See Documentation: The radius should be set to at least the size of the largest object that is not part of the background.ss
	saveAs("Tiff", filenameNew);}
	close();}
	