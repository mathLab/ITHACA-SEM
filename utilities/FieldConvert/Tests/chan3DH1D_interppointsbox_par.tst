<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interp field to a box of points (also calculate cp and cp0)</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e  -m interppoints:cp=0,0.5:box=10,10,9,0,1,0,1,0,1:fromxml=chan3DH1D.xml:fromfld=chan3DH1D_zVariation.fld chan3DH1D_box.dat</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
        <file description="Session File">chan3DH1D_zVariation.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6"> 0.593171</value>
            <value variable="y" tolerance="1e-6">0.593171</value>
            <value variable="z" tolerance="1e-6">0.493007</value>
            <value variable="u" tolerance="1e-6">1.73509</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">0</value>
            <value variable="Cp" tolerance="1e-6">0</value>
            <value variable="Cp0" tolerance="1e-6">3.58535</value>
        </metric>
    </metrics>
</test>