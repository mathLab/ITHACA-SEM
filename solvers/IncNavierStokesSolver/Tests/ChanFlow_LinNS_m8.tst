<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Unsteady channel flow with coupled solver , P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_LinNS_m8.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_LinNS_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">4.07942e-14</value>
            <value variable="v" tolerance="1e-6">2.91568e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">7.62834e-13</value>
            <value variable="v" tolerance="1e-6">2.14091e-13</value>
        </metric>
    </metrics>
</test>


