package users.veronika;

import de.embl.cba.morphometry.algorithm.DirectionalLineFilter;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij3d.ContentConstants;
import ij3d.Image3DUniverse;
import inra.ijpb.morphology.directional.DirectionalFilter;

public class ExploreDirectionalLineFilter
{
	public static void main( String[] args )
	{
		final ImageJ imageJ = new ImageJ();

		final ImagePlus imagePlus = IJ.openImage( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/users/veronika/straight-lines-test.zip" );

		final DirectionalLineFilter filter = new DirectionalLineFilter();
		filter.run( imagePlus, DirectionalFilter.Operation.MEAN, 60, 20 );
		final ImagePlus result = filter.getResult();
		result.show();
		result.getCalibration().pixelDepth = 10;

		final Image3DUniverse universe = new Image3DUniverse();
		universe.show();
		universe.addContent( result, ContentConstants.VOLUME );

	}
}
