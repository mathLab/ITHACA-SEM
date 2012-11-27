<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.15938e-16</value>
            <value variable="v" tolerance="1e-12">1.73924e-16</value>
	    <value variable="p" tolerance="1e-12">9.13284e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.50435e-14</value>
            <value variable="v" tolerance="1e-12">1.0547e-15</value>
	    <value variable="p" tolerance="1e-12">8.199e-14</value>
        </metric>
    </metrics>
</test>


