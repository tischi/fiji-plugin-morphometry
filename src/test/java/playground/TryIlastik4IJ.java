package playground;

import bdv.viewer.Source;
import bdv.viewer.SourceAndConverter;
import mpicbg.spim.data.SpimData;
import org.ilastik.ilastik4ij.IlastikOptions;
import org.ilastik.ilastik4ij.IlastikPixelClassificationPrediction;

public class TryIlastik4IJ
{
	public static void main( String[] args )
	{
		final IlastikPixelClassificationPrediction prediction = new IlastikPixelClassificationPrediction();
		final IlastikOptions options = new IlastikOptions();
		options.setExecutableFilePath( "/Applications/" );

//		new SpimData(  )

	}
}
