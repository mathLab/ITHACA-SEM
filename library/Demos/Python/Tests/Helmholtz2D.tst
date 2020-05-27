<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz solver in 2D domain</description>
    <executable python="true"> Helmholtz2D.py </executable>
    <parameters>Helmholtz2D_P7.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">../../MultiRegions/Tests/Helmholtz2D_P7.xml</file>
    </files>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-7" variable="nek">6.120209e-05</value>
            <value tolerance="1e-7" variable="nekpy">6.120209e-05</value>
        </metric>
        <metric type="L2" id="2">
            <value tolerance="1e-7" variable="nek">4.302677e-05</value>
        </metric>
        <metric type="regex" id="3">
            <regex>^Reduction test.*: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">4</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
