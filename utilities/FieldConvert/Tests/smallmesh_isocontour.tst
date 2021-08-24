<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Extract an isocontour</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -n6 -m isocontour:fieldstr="p":fieldvalue=0.1:globalcondense:smooth smallmesh.xml smallmesh.fld iso.dat</parameters>
    <files>
        <file description="Session File">smallmesh.xml</file>
        <file description="Session File">smallmesh.fld</file>
    </files>

     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.644241</value>
            <value variable="y" tolerance="1e-4">0.450692</value>
            <value variable="z" tolerance="1e-4">0.434032</value>
            <value variable="u" tolerance="1e-4">0.901556</value>
            <value variable="v" tolerance="1e-4">0.101557</value>
            <value variable="w" tolerance="1e-4">0.0892263</value>
            <value variable="p" tolerance="1e-4">0.1</value>
            <value variable="isocon" tolerance="1e-4">0.1</value>
        </metric>
    </metrics>
</test>

