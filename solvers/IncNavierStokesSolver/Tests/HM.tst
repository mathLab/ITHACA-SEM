<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fourier Half Mode Basis, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>HM.xml</parameters>
    <files>
        <file description="Session File">HM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.48637e-16</value>
            <value variable="v" tolerance="1e-12">1.90061e-16</value>
            <value variable="w" tolerance="1e-12">1.566e-16</value>
            <value variable="p" tolerance="1e-12">1.27194e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.11022e-15</value>
            <value variable="v" tolerance="1e-12">9.71445e-16</value>
            <value variable="w" tolerance="1e-12">6.10623e-16</value>
            <value variable="p" tolerance="1e-12">4.89839e-14</value>
        </metric>
    </metrics>
</test>


