<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interplation on a sphere surface </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=interpfield_fromsphere.xml:fromfld=interpfield_fromsphere.fld interpfield_tosphere.xml tosphere.dat</parameters>
    <files>
        <file description="Tri and quad elements">interpfield_fromsphere.xml</file>
        <file description="Tri elements">interpfield_tosphere.xml</file>
        <file description="Field File">interpfield_fromsphere.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">2.04663</value>
            <value variable="y" tolerance="1e-4">2.04662</value>
            <value variable="z" tolerance="1e-4">2.04662</value>
            <value variable="eta" tolerance="1e-4">1.15612</value>
            <value variable="u" tolerance="1e-4">2.04663</value>
            <value variable="v" tolerance="1e-4">3.54488</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">1.</value>
            <value variable="y" tolerance="1e-4">1.</value>
            <value variable="z" tolerance="1e-4">1.</value>
            <value variable="eta" tolerance="1e-4">0.770749</value>
            <value variable="u" tolerance="1e-4">1.</value>
            <value variable="v" tolerance="1e-4">1.</value>
        </metric>
    </metrics>
</test>
