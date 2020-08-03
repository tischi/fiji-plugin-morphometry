package microglia;

import de.embl.cba.morphometry.microglia.MicrogliaMorphometryCommand;
import net.imagej.ImageJ;

import java.io.File;

public class TestMicrogliaMorphometry
{
	//@Test
	public void test()
	{
		final ImageJ imageJ = new ImageJ();
		final MicrogliaMorphometryCommand command = new MicrogliaMorphometryCommand();
		command.showIntermediateResults = false;
		command.opService = imageJ.op();
		command.outputDirectory = new File("/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/test-data/microglia");
		command.labelMaskFile = new File("/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/test-data/microglia/MAX_5C-crop-t1-3-intensities-labelMasks.tif");
		command.intensityFile = new File("/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/test-data/microglia/MAX_5C-crop-t1-3-intensities.tif");
		command.run();
	}

	public static void main( String[] args )
	{
 		new TestMicrogliaMorphometry().test();
	}
}
