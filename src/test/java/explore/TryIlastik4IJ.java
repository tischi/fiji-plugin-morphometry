package explore;

import bdv.viewer.Source;
import bdv.viewer.SourceAndConverter;
import ij.IJ;
import ij.ImagePlus;
import mpicbg.spim.data.SpimData;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import org.ilastik.ilastik4ij.IlastikOptions;
import org.ilastik.ilastik4ij.IlastikPixelClassificationPrediction;

import java.io.File;

public class TryIlastik4IJ
{
	public static < R extends RealType< R > > void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		final IlastikOptions options = new IlastikOptions();
		options.setExecutableFilePath( "/Applications/ilastik-1.3.3-OSX.app" );
		options.setNumThreads( 4 );
		options.setMaxRamMb( 10000 );

		final ImagePlus imagePlus = IJ.openImage( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/test-data/spindle/SpindleVolumeTest01.zip" );

		final RandomAccessibleInterval< R > raiXYCZ = ImageJFunctions.wrapReal( imagePlus );

		DatasetService datasetService = ij.dataset();
		AxisType[] axisTypes = new AxisType[]{ Axes.X, Axes.Y, Axes.CHANNEL, Axes.Z };
		ImgPlus imgPlus = new ImgPlus( datasetService.create( raiXYCZ  ), "image", axisTypes );
		ij.ui().show( imgPlus );


		final IlastikPixelClassificationPrediction prediction = new IlastikPixelClassificationPrediction();
		prediction.log = ij.log();
		prediction.statusService = ij.status();
		prediction.optionsService = ij.options();
		prediction.chosenOutputType = "Segmentation";
		prediction.ilastikOptions = options;
		prediction.inputImage = datasetService.create( imgPlus );
		prediction.projectFileName = new File("/Users/tischer/Documents/tobias-kletter/20191206_DNA_Segmentation_2Ch.ilp");
		prediction.run();
	}
}
