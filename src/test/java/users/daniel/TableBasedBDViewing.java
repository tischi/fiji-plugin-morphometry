package users.daniel;

import de.embl.cba.morphometry.commands.BDOpenTableCommand;
import net.imagej.ImageJ;

import javax.swing.*;
import java.io.File;

public class TableBasedBDViewing
{
	public static void main( String[] args )
	{
		final ImageJ imageJ = new ImageJ();
		imageJ.ui().showUI();
		final BDOpenTableCommand command = new BDOpenTableCommand();
		command.imageTablePath = new File("/Users/tischer/Documents/BD-image-processing/sample_data/minimalgated/countpath.csv");
		command.commandService = imageJ.command();
		command.run();
	}
}
