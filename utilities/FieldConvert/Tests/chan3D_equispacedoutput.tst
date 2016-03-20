<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a field into a equi-spaced tecplot file</description>
    <executable>FieldConvert</executable>
    <parameters>-m equispacedoutput -e chan3D.xml chan3D.fld equispacedoutput.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.745356</value>
            <value variable="y" tolerance="1e-4">0.745356</value>
            <value variable="z" tolerance="1e-4">0.745356</value>
            <value variable="u" tolerance="1e-4">0.628539</value>
            <value variable="v" tolerance="1e-4">0</value>
            <value variable="w" tolerance="1e-4">0</value>
            <value variable="p" tolerance="1e-4">2.49444</value>
        </metric>
    </metrics>
</test>

