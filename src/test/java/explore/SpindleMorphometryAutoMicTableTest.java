package explore;//import automic.table.TableModel;
//import automic.utils.imagefiles.PatternOpenerWithBioformats;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

public class SpindleMorphometryAutoMicTableTest<T extends RealType<T> & NativeType< T > >
{

	static final String highZoomPatternTag = "High.Zoom";

	public void run() throws Exception
	{
//		//specify folder with your ImageJ plugins. This allows IJ.run("Bio-Formats",...) to work
//		System.getProperties().setProperty("plugins.dir", "/Applications/Fiji.app/plugins");
//
//		ImageJ imagej = new ImageJ();
//		imagej.ui().showUI();
//
//		//create table model from labelMaskFile
//		TableModel tbl = new TableModel( new File("/Volumes/almfscreen/kletter/20180904/formattedLog_0906-1_Sorted_Cropped.txt") );
//
//		//get number of datasets in the table
//		int nDatasets = tbl.getRowCount();
//
//		//loop through datasets
//		for ( int iDataset=0; iDataset<nDatasets; iDataset++ )
//		{
//			if ( !tbl.getBooleanValue( iDataset, "Success" ) )//no high zoom image acquired
//				continue;
//
//			//get folder and filename components of the pattern from the table
//			String patternPath = tbl.getFileParentPathString( iDataset, highZoomPatternTag, "IMGR" );
//			String patternName = tbl.getFileName( iDataset, highZoomPatternTag, "IMGR" );
//
//			final String absolutePath = new File( patternPath, patternName ).getAbsolutePath();
//
//			ImagePlus imp = PatternOpenerWithBioformats.openImage( absolutePath );
//
//			final RandomAccessibleInterval< T > wrap = ImageJFunctions.wrap( imp );
//			final RandomAccessibleInterval< T > dnaImage = Views.hyperSlice( wrap, 2, 1);
//			final RandomAccessibleInterval< T > tubulinImage = Views.hyperSlice( wrap, 2, 0);
//
//			SpindleMorphometrySettings settings = new SpindleMorphometrySettings();
//			settings.dnaImage = dnaImage;
//			settings.tubulinImage = tubulinImage;
//			settings.showIntermediateResults = false;
//			settings.inputCalibration = Utils.getCalibration( imp );
//			settings.inputCalibration[ 2 ] = 1.05; // TODO: UI
//			settings.workingVoxelSize = 0.25;
//			settings.maxPossibleValueInDataSet = Math.pow( 2, imp.getBitDepth() ) - 1.0;
//			settings.maxShortAxisDist = 6;
//			settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
//			settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
//			settings.watershedSeedsGlobalDistanceThreshold = 2.0;
//			settings.interestPointsRadius = 0.25;
//			String outDir = "/Volumes/almfscreen/kletter/20180904/tischi";
//			final String[] split = new File( patternPath ).toString().split( "/" );
//			final String subfolder = split[ split.length - 4 ] + "--" + split[ split.length - 2 ] + "--" + split[ split.length - 1 ];
//			outDir += File.separator + subfolder;
//			settings.outputDirectory = new File( outDir );
//			settings.inputDataSetName = patternName;
//
//			Utils.log( "Analyzing: " + subfolder + "/" + patternName );
//			Utils.log( "Output directory: " + settings.outputDirectory );
//
//			SpindleMorphometry morphometry = new SpindleMorphometry( settings, imagej.op() );
//			morphometry.run();
//		}

	}

	public static void main( String... args )
	{
		try
		{
			new SpindleMorphometryAutoMicTableTest().run();
		}
		catch ( Exception e )
		{
			e.printStackTrace();
		}

	}
}
