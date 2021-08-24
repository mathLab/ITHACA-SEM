<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Project3D Pyramid Modified basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 6 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.5 0.5 0.866</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">1.15445e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">8.91376e-12</value>
        </metric>
    </metrics>
</test>
