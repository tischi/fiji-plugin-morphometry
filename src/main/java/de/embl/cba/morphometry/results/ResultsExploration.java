package de.embl.cba.morphometry.results;

import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.tables.InteractiveTablePanel;
import ij.IJ;
import ij.ImagePlus;
import org.scijava.table.GenericTable;

import java.io.File;
import java.util.ArrayList;

public class ResultsExploration
{

	private void showResults( File file, ArrayList< String > measurements )
	{
		final ImagePlus imagePlus = IJ.openImage( file.getAbsolutePath() );
		imagePlus.show();

//		final GenericTable table = Measurements.createGenericTableFromTableRows( measurements );
//
//		int[] xyzt = new int[ 4 ];
//		xyzt[ 0 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + "X" + Measurements.SEP + Measurements.PIXEL_UNITS );
//		xyzt[ 1 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + "Y" + Measurements.SEP + Measurements.PIXEL_UNITS );
//		xyzt[ 2 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + "Z" + Measurements.SEP + Measurements.PIXEL_UNITS );
//		xyzt[ 3 ] = table.getColumnIndex( Measurements.COORDINATE + Measurements.SEP + Measurements.TIME + Measurements.SEP +  Measurements.FRAME_UNITS );
//
//		final InteractiveTablePanel interactiveTablePanel = new InteractiveTablePanel( table );
//		interactiveTablePanel.setCoordinateColumns( xyzt );
//		interactiveTablePanel.setImagePlus( imagePlus );
	}
}
