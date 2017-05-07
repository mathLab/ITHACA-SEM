<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interp field to a box of points (also calculate cp and cp0)</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e  -m interppoints:cp=0,0.5:box=10,10,10,-0.5,0.5,-0.5,0.5,-0.5,0.5:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_box.dat</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.319142</value>
            <value variable="y" tolerance="1e-6">0.319142</value>
            <value variable="z" tolerance="1e-6">0.319142</value>
            <value variable="u" tolerance="1e-6">0.902617</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">2.09938</value>
            <value variable="Cp" tolerance="1e-6">4.19877</value>
            <value variable="Cp0" tolerance="1e-6">4.98355</value>
        </metric>
    </metrics>
</test>

