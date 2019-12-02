package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.fccf.FCCF;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import loci.common.DebugTools;
import org.scijava.command.Command;
import org.scijava.command.CommandService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.widget.Button;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>FCCF>BD Processing" )
public class BDOpenTableCommand implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public CommandService commandService;

	@Parameter ( label = "Image Table" )
	public File imageTablePath;

	public void run()
	{
		// load table
		final ArrayList< String > classes = new ArrayList<>();
		classes.add( "Hello" );
		classes.add( "World" );
		BDImageViewingCommand.classChoices = classes;
		commandService.run( BDImageViewingCommand.class, true );
	}


}
