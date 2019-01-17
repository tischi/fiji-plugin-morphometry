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

public class TranslocationTestCommand
{
	public static  < T extends RealType< T > & NativeType< T > >
	void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		ij.ui().showUI();

		final ImagePlus imagePlus = IJ.openImage(
				TranslocationTestCommand.class.getResource(
						"translocation/test01.zip" ).getFile() );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open",
				TranslocationTestCommand.class.getResource(
						"translocation/test01-rois.zip" ).getFile() );
		final ArrayList< FinalInterval > intervals = Rois.asIntervals( rm.getRoisAsArray() );

		ArrayList< RandomAccessibleInterval< T > > intensities = Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 );

		final MembraneTranslocationComputer computer = new MembraneTranslocationComputer(
				intensities,
				intervals,
				ij.op() );

		final ArrayList< TranslocationResult > results = computer.getResults();

		for ( int r = 0; r < results.size(); r++ )
		{
			Utils.listOf2DImagesAsImagePlusMovie( results.get( r ).binaryGradients, ""+(r+1)+"-binary-gradients" ).show();
			Utils.listOf2DImagesAsImagePlusMovie( results.get( r ).intensities, ""+(r+1)+"-intensities" ).show();
			Utils.listOf2DImagesAsImagePlusMovie( results.get( r ).gradients, ""+(r+1)+"-gradients" ).show();
			Utils.listOf2DImagesAsImagePlusMovie( results.get( r ).membranes, ""+(r+1)+"-membranes" ).show();
		}

		final JTable table = TranslocationResult.resultsAsTable( results );

		final ObjectTablePanel objectTablePanel = new ObjectTablePanel( table );

		objectTablePanel.showPanel();

		final ArrayList< RandomAccessibleInterval< T > > labelMasks =
				TranslocationCommand.createLabelMasks( intensities, results );

		TranslocationCommand.showLabelMasksAndIntensities( intensities, labelMasks, imagePlus.getTitle() );

		TranslocationCommand.plotTranslocations( results );

	}
}
