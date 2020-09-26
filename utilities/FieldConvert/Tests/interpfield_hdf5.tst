<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interp field to a box of points (also calculate cp and cp0) for HDF5 input/output</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=interptest_hdf5.xml:fromfld=interptest_hdf5.fld  interptest.xml out.fld </parameters>
    <files>
        <file description="Session File">interptest.xml</file>
        <file description="Session File">interptest_hdf5.xml</file>
        <file description="Session File">interptest_hdf5.nekg</file>
        <file description="Session File">interptest_hdf5.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.0388428</value>
            <value variable="y" tolerance="1e-4">0.0157524</value>
            <value variable="u" tolerance="1e-4">0.0254815</value>
            <value variable="v" tolerance="1e-4">0.00331392</value>
            <value variable="p" tolerance="1e-4">0.00229271</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">0.967848</value>
            <value variable="y" tolerance="1e-4">0.407481</value>
            <value variable="u" tolerance="1e-4">0.784991</value>
            <value variable="v" tolerance="1e-4">0.390773</value>
            <value variable="p" tolerance="1e-4">0.0871031</value>
        </metric>
    </metrics>
</test>

