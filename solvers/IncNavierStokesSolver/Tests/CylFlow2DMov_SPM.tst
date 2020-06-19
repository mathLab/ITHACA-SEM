<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D cylinder flow, P=3, cylinder defined via IB SPM function</description>    
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow2DMov_SPM.xml</parameters>
    <files>
        <file description="Session File">CylFlow2DMov_SPM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">22.3913</value>
            <value variable="v" tolerance="1e-6">0.707402</value>
            <value variable="p" tolerance="1e-6">7.23147</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">2.02754</value>
            <value variable="v" tolerance="1e-6">0.94615</value>
            <value variable="p" tolerance="1e-6">4.16306</value>
        </metric>
    </metrics>
</test>
