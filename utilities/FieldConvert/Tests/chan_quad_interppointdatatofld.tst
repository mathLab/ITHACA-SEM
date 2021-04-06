<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interpolate a .csv onto a mesh </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interppointdatatofld:frompts=chan_quad_interppointdatatofld.csv chan_quad_interppointdatatofld.xml chan_quad_interppointdatatofld.fld</parameters>
    <files>
        <file description="Session File">chan_quad_interppointdatatofld.xml</file>
        <file description="Session File">chan_quad_interppointdatatofld.csv</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">2.12132</value>
            <value variable="y" tolerance="1e-4">1.41421</value>
            <value variable="u" tolerance="1e-3">3.22044</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">1.5</value>
            <value variable="y" tolerance="1e-4">1</value>
            <value variable="u" tolerance="1e-3">3.23248</value>
        </metric>
    </metrics>
</test>
