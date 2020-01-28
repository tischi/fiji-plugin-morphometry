package de.embl.cba.morphometry.commands;

import de.embl.cba.tables.Tables;
import ij.IJ;
import org.scijava.command.Command;
import org.scijava.command.CommandService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.util.List;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>FCCF>BD>Process Images from Folder" )
public class BDOpenFolderCommand implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public CommandService commandService;

	@Parameter ( label = "Input images directory" )
	public File inputImagesDirectory;

	public void run()
	{
		BDImageViewingAndSavingCommand.jTable = null;
		BDImageViewingAndSavingCommand.inputImagesDirectory = inputImagesDirectory;

		commandService.run( BDImageViewingAndSavingCommand.class, true );
	}
}
