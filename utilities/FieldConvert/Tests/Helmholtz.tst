<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process gradient </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m gradient Helmholtz.xml Helmholtz.fld Helmholtz.dat</parameters>
    <files>
        <file description="Session File">Helmholtz.xml</file>
        <file description="Session File">Helmholtz.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x"   tolerance="1e-4">1.1547</value>
            <value variable="y"   tolerance="1e-4">1.1547</value>
            <value variable="u"   tolerance="1e-4">1.00002</value>
            <value variable="u_x" tolerance="1e-4">3.14172</value>
            <value variable="u_y" tolerance="1e-4">3.14163</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"   tolerance="1e-4">1.</value>
            <value variable="y"   tolerance="1e-4">1.</value>
            <value variable="u"   tolerance="1e-4">1.</value>
            <value variable="u_x" tolerance="1e-4">3.14199</value>
            <value variable="u_y" tolerance="1e-4">3.12754</value>
        </metric>
    </metrics>
</test>

