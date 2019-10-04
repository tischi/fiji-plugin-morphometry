run("Close All");

targetVoxelSize = 0.25 // um

inputDir = "/Users/tischer/Desktop/3d_Plugin_Problematic_Images";

files = newArray(1);
files[ 0 ] = "R1ENestin(h)mCherry_doubletag_Unfacs_Tubb3KO_clone16_046_Cropped_1_ExceptionDuringComputation.tif";
//files[ 0 ] = "aaa.tif";

for (i = 0; i < files.length; i++) 
{
	path = inputDir + "/" + files[ i ];
	print( "Opening: " + path );
	open( path ); 
	rename( "image" );
	scale( "image" );
}

function scale( name ) {
	selectWindow( name );
	
	getVoxelSize(width, height, depth, unit);
	scaleX = 1.0 * width / targetVoxelSize;
	scaleY = 1.0 * height / targetVoxelSize;
	scaleZ = 1.0 * depth / targetVoxelSize;
	getDimensions(width, height, channels, slices, frames);
	newSlices = slices * scaleZ;
	print( scaleZ );
	run("Scale...", "x=&scaleX y=&scaleY z=&scaleZ depth=&newSlices interpolation=Bilinear average create");
	rename( "scaled" );
	selectWindow( name ); close( );
	selectWindow( "scaled" );
	rename( name );
}
