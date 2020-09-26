<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D cylinder flow, P=3, cylinder defined via IB SPM function</description>    
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow2D_SPM.xml</parameters>
    <files>
        <file description="Session File">CylFlow2D_SPM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">22.4052</value>
            <value variable="v" tolerance="1e-6">0.70658</value>
            <value variable="p" tolerance="1e-6">7.64851</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1.77878</value>
            <value variable="v" tolerance="1e-6">1.06215</value>
            <value variable="p" tolerance="1e-6">3.57015</value>
        </metric>
    </metrics>
</test>
