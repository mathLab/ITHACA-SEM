<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project3D Prism Modified basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.5 0.0 1.0 0.5 1.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">5.16422e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">2.08524e-11</value>
        </metric>
    </metrics>
</test>
