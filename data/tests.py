from django.test import TestCase


class StaticPagesSmokeTest(TestCase):
    """
    Basic smoke tests: these pages are plain TemplateViews with no
    Galaxy/DB dependency, so they must render successfully in any
    environment (including CI, which has no Galaxy server or Redis).
    """

    def test_pages_return_200(self):
        paths = [
            '/',
            '/about',
            '/documentation',
            '/analysis',
            '/status',
            # Exercises real crispy_forms + django-simple-captcha
            # rendering (not just app loading), unlike the other pages.
            '/about/feedback',
        ]
        for path in paths:
            response = self.client.get(path)
            self.assertEqual(
                response.status_code, 200,
                "GET %s returned %s" % (path, response.status_code))
