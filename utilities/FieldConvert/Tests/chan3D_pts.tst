<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .fld </description>
    <executable>FieldConvert</executable>
    <parameters>-e -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_pts.csv out.csv</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
        <file description="Session File">chan3D_pts.csv</file>
    </files>
    <!-- currently no metrics supported in OutputPts, we only test if it crashes -->
</test>

