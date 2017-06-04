# perl_modules

Modules with commonly used subroutines

To use any of these subroutines in a Perl script, add a module file (e.g., sloan.pm) to a location that is in the set of paths that Perl searches for modules. Locations can be added by modifying the PERL5LIB environmental variable.

For example, in a Mac or Linux system, add the following line to a .bashrc or .bash_profile file in your home directory:

export PERL5LIB=/PATH/TO/NEW/DIR:$PERL5LIB

(where /PATH/TO/NEW/DIR is the actual path to the directory with the module, not all the way to the file itself)

Then simply add the line "use sloan;" to your Perl script, and all the functions in that module will be available for use in your script.
