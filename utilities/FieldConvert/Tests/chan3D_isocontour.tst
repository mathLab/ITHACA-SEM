<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Extract an isocontour</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m isocontour:fieldstr="u+v":fieldvalue=0.5:fieldname="UplusV":globalcondense:smooth chan3D.xml chan3D.fld isocontour.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>

     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.597226</value>
            <value variable="y" tolerance="1e-4">0.595038</value>
            <value variable="z" tolerance="1e-4">0.666667</value>
            <value variable="u" tolerance="1e-4">0.5</value>
            <value variable="v" tolerance="1e-4">0</value>
            <value variable="w" tolerance="1e-4">0</value>
            <value variable="p" tolerance="1e-4">2.36834</value>
            <value variable="UplusV" tolerance="1e-4">0.5</value>
        </metric>
    </metrics>
</test>

