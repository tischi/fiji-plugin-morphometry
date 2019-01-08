import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.commands.TranslocationCommand;
import de.embl.cba.morphometry.translocation.MembraneTranslocationComputer;
import de.embl.cba.morphometry.translocation.TranslocationResult;
import de.embl.cba.tables.objects.ObjectTablePanel;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import javax.swing.*;
import java.util.ArrayList;

public class TranslocationExample
{
	public static  < T extends RealType< T > & NativeType< T > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		ij.ui().showUI();

		final ImagePlus imagePlus = IJ.openImage(
				TranslocationExample.class.getResource(
						"translocation/test03.zip" ).getFile() );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open",
				TranslocationExample.class.getResource(
						"translocation/test03.roi" ).getFile() );
		final ArrayList< FinalInterval > intervals = Rois.asIntervals( rm.getRoisAsArray() );

		ArrayList< RandomAccessibleInterval< T > > intensities = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MembraneTranslocationComputer computer = new MembraneTranslocationComputer(
				intensities,
				intervals,
				ij.op() );

		final ArrayList< TranslocationResult > results = computer.getResults();

		for ( int r = 0; r < results.size(); r++ )
		{
			Utils.listOf2DImagesAsImagePlusMovie( results.get( r ).gradients, "gradients" ).show();
			Utils.listOf2DImagesAsImagePlusMovie( results.get( r ).membraneMasks, "membraneMasks" ).show();
		}

		final JTable table = TranslocationResult.resultsAsTable( results );

		final ObjectTablePanel objectTablePanel = new ObjectTablePanel( table );

		objectTablePanel.showPanel();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks =
				TranslocationCommand.createLabelMasks( intensities, results );

		TranslocationCommand.showLabelMasksAndIntensities( intensities, labelMasks );

		TranslocationCommand.plotTranslocations( results );

	}
}
