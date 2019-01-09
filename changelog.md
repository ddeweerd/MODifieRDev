# Changelog

### Version 0.1.2

`added` Settings list to input objects

`added` Function `recalculate_diff_genes` to update existing input object DEGs
to change to/from adjusted p values.

`added` Function `recalculate_expression` to update existing input_object to
change the method of collapsing probe_names

`changed` The settings list in the objects to evaluate specific arguments  
For example, if a variable has been used as a parameter, the variable will be
evaluated and the value will appear in the settings instead of the variable name.

`fixed` The annotation table was missing from the input object, has been fixed

`fixed` Module set object gene_by_method counting
