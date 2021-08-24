<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>C0 Projection 2D with P=7</description>
    <executable>ProjectCont2D</executable>
    <parameters>ProjectionCont2D_P7.xml</parameters>
    <files>
        <file description="Session File">ProjectionCont2D_P7.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">2.10047e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">2.28465e-09</value>
        </metric>
    </metrics>
</test>


